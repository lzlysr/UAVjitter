// 独立UPA积分诊断的可选MEX加速内核，不改变节点、Capon公式或求积权重。
// 仅将MATLAB分块中的有限傅里叶求和改为本地并行循环；所有累计使用双精度。
// 输入：系数、阵元差、两个坐标轴、规则编号、阵列参数、名义方向、排除半宽、
// 测度编号和线程数。输出：三角矩、纳入节点数、积分域测度。
#include "mex.h"
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <atomic>
#include <omp.h>

extern "C" __declspec(dllexport) void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]) {
    if(nrhs!=10 || nlhs!=3)
        mexErrMsgIdAndTxt("uav:incm:NativeArguments","本地诊断内核需要10个输入和3个输出。");
    const double* coeff=mxGetPr(prhs[0]);
    const int K=static_cast<int>(mxGetM(prhs[0])),B=static_cast<int>(mxGetN(prhs[0]));
    const int Q=static_cast<int>(mxGetM(prhs[1]));
    const double* lags=mxGetPr(prhs[1]);
    const double* axis1=mxGetPr(prhs[2]);const double* axis2=mxGetPr(prhs[3]);
    const int G=static_cast<int>(mxGetNumberOfElements(prhs[2]));
    const int rule=static_cast<int>(mxGetScalar(prhs[4]));
    const double* array=mxGetPr(prhs[5]);const int Nx=static_cast<int>(array[0]),Ny=static_cast<int>(array[1]);
    const double* reference=mxGetPr(prhs[6]);const double beta=mxGetScalar(prhs[7]);
    const double* measures=mxGetPr(prhs[8]);const int C=static_cast<int>(mxGetNumberOfElements(prhs[8]));
    const int threads=static_cast<int>(mxGetScalar(prhs[9]));
    if(K!=1+2*Q || mxGetN(prhs[1])!=2 || G<3 || B<1 || C<1 || threads<1 ||
       Nx<2 || Ny<2 || mxGetNumberOfElements(prhs[3])!=static_cast<mwSize>(G))
        mexErrMsgIdAndTxt("uav:incm:NativeShape","本地诊断内核输入尺寸无效。");
    const double pi=std::acos(-1.0),threshold=std::cos(beta);
    // 系数转为按试次连续存储，便于编译器对小批量协方差执行向量化。
    std::vector<double> cf(K*B);
    for(int k=0;k<K;k++)for(int b=0;b<B;b++)cf[k*B+b]=coeff[k+K*b];
    std::vector<std::vector<double>> partial(threads,std::vector<double>(K*B*C,0.0));
    std::vector<std::vector<double>> mass(threads,std::vector<double>(C,0.0));
    std::vector<unsigned long long> counts(threads,0);
    std::atomic<bool> invalid(false);
    const long long total=static_cast<long long>(G)*G;
    #pragma omp parallel num_threads(threads)
    {
        const int tid=omp_get_thread_num();
        std::vector<std::complex<double>> px(Nx),py(Ny);
        std::vector<double> basis(K),den(B),power(B),weights(C);
        double* acc=partial[tid].data();
        #pragma omp for schedule(static)
        for(long long n=0;n<total;n++) {
            const int row=static_cast<int>(n%G),col=static_cast<int>(n/G);
            double x,y,z=0,r=0;
            if(rule==0 || rule==1) {x=axis1[col];y=axis2[row];}
            else {z=axis2[row];r=std::sqrt(std::max(0.0,1-z*z));x=r*std::cos(axis1[col]);y=r*std::sin(axis1[col]);}
            const double rr=x*x+y*y;
            if(rr>1)continue;
            const double restoredZ=-std::sqrt(std::max(0.0,1-rr));
            if(x*reference[0]+y*reference[1]+restoredZ*reference[2]>=threshold)continue;
            counts[tid]++;
            for(int c=0;c<C;c++) {
                if(rule==0)weights[c]=4.0/((G-1.0)*(G-1.0));
                else if(rule==1)weights[c]=4.0/(static_cast<double>(G)*G);
                else if(rule==2)weights[c]=2*pi/(static_cast<double>(G)*G)*(measures[c]==1?std::abs(z):1.0);
                else weights[c]=pi*pi/(static_cast<double>(G)*G)*(measures[c]==3?1.0:r*(measures[c]==1?std::abs(z):1.0));
                mass[tid][c]+=weights[c];
            }
            const double ax=2*pi*array[2]*x,ay=2*pi*array[3]*y;
            const std::complex<double> ex(std::cos(ax),std::sin(ax)),ey(std::cos(ay),std::sin(ay));
            px[0]=1.0;py[0]=1.0;
            for(int k=1;k<Nx;k++)px[k]=px[k-1]*ex;
            for(int k=1;k<Ny;k++)py[k]=py[k-1]*ey;
            basis[0]=1.0;
            for(int q=0;q<Q;q++) {
                const int lx=static_cast<int>(lags[q]),ly=static_cast<int>(lags[q+Q]);
                const auto value=px[lx]*(ly<0?std::conj(py[-ly]):py[ly]);
                basis[1+q]=value.real();basis[1+Q+q]=value.imag();
            }
            std::fill(den.begin(),den.end(),0.0);
            for(int k=0;k<K;k++) {
                const double v=basis[k];const double* ck=cf.data()+k*B;
                #pragma omp simd
                for(int b=0;b<B;b++)den[b]+=v*ck[b];
            }
            bool ok=true;for(int b=0;b<B;b++)if(!(den[b]>0) || !std::isfinite(den[b]))ok=false;
            if(!ok){invalid=true;continue;}
            for(int c=0;c<C;c++) {
                for(int b=0;b<B;b++)power[b]=weights[c]/den[b];
                for(int k=0;k<K;k++) {
                    double* out=acc+(c*K+k)*B;const double v=basis[k];
                    #pragma omp simd
                    for(int b=0;b<B;b++)out[b]+=v*power[b];
                }
            }
        }
    }
    if(invalid)mexErrMsgIdAndTxt("uav:incm:NativeDenominator","本地内核遇到非正或非有限Capon分母。");
    const mwSize dims[3]={static_cast<mwSize>(K),static_cast<mwSize>(B),static_cast<mwSize>(C)};
    plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1]=mxCreateDoubleScalar(0);plhs[2]=mxCreateDoubleMatrix(1,C,mxREAL);
    double* out=mxGetPr(plhs[0]);double* outmass=mxGetPr(plhs[2]);
    unsigned long long count=0;
    for(int t=0;t<threads;t++) {
        count+=counts[t];
        for(int c=0;c<C;c++) {
            outmass[c]+=mass[t][c];
            for(int k=0;k<K;k++)for(int b=0;b<B;b++)out[k+K*b+K*B*c]+=partial[t][(c*K+k)*B+b];
        }
    }
    *mxGetPr(plhs[1])=static_cast<double>(count);
}
