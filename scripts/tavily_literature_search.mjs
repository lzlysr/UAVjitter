// 调用本机已配置版本的 Tavily MCP，保存原始结果以供文献筛选复核。
// 不记录或输出密钥；检索结果只作为发现线索，公式与阵列类型仍须查正文。
import {createRequire} from 'node:module';
import {pathToFileURL} from 'node:url';
import {mkdir, writeFile} from 'node:fs/promises';
import path from 'node:path';
import {createInterface} from 'node:readline';

const packageRoot = process.env.TAVILY_PACKAGE_ROOT;
if (!packageRoot) throw new Error('请用 TAVILY_PACKAGE_ROOT 指向已安装的 tavily-mcp 目录。');
const require = createRequire(path.join(packageRoot, 'package.json'));
const {Client} = await import(pathToFileURL(require.resolve('@modelcontextprotocol/sdk/client/index.js')));
const {StdioClientTransport} = await import(pathToFileURL(require.resolve('@modelcontextprotocol/sdk/client/stdio.js')));
const client = new Client({name:'upa-ieee-literature-audit',version:'1.0.0'});
// MCP SDK 默认只传递一组安全环境变量；必须显式传入会话中的认证变量。
const transport = new StdioClientTransport({command:process.execPath,args:[path.join(packageRoot,'build/index.js')],env:process.env,stderr:'pipe'});
const [mode, ...arguments_] = process.argv.slice(2);
try {
  await client.connect(transport);
  if (mode === 'session') {
    console.log('Tavily 会话就绪；已启用密钥认证:',Boolean(process.env.TAVILY_API_KEY));
    for await (const line of createInterface({input:process.stdin})) {
      if (line.trim() === 'exit') break;
      try {
        const input = JSON.parse(line);
        await run(input.mode ?? 'batch', input.queries ?? []);
      } catch (error) {
        const message = String(error.message);
        console.log('请求失败:',process.env.TAVILY_API_KEY ? message.replaceAll(process.env.TAVILY_API_KEY,'[已隐藏]') : message);
      }
      console.log('TAVILY_REQUEST_COMPLETE');
    }
  } else await run(mode, arguments_);
} finally {
  await client.close();
  // Windows 交互终端的 stdin 可能继续保持事件循环，关闭后不再挂住认证会话。
  if (mode === 'session') {
    process.stdin.pause();
    process.stdin.unref?.();
  }
}

async function run(mode, arguments_) {
  if (mode === 'list') {
    console.log(JSON.stringify(await client.listTools(),null,2));
  } else {
    const queries = mode === 'batch' ? arguments_ : [arguments_.join(' ')];
    for (const query of queries) {
    const response = await client.callTool({name:mode === 'extract' ? 'tavily_extract' : 'tavily_search',arguments:
      mode === 'extract' ? {urls:query.split(' '),extract_depth:'advanced',format:'text'} :
      {query,search_depth:'advanced',max_results:10,include_raw_content:true,include_answer:false}},undefined,{timeout:90000});
    // 防止第三方错误信息意外回显凭据；保存和展示前再次脱敏。
    const serialized = JSON.stringify(response);
    const result = JSON.parse(process.env.TAVILY_API_KEY ? serialized.replaceAll(process.env.TAVILY_API_KEY,'[已隐藏]') : serialized);
    const out = path.resolve('data/literature/tavily_upa_ieee_2026-10-20');
    await mkdir(out,{recursive:true});
    const filename = path.join(out,`${Date.now()}_${mode}.json`);
    await writeFile(filename,JSON.stringify({mode,query,result},null,2));
    console.log('原始记录:',filename);
    const content = result.content?.filter(x=>x.type==='text').map(x=>x.text).join('\n') ?? JSON.stringify(result);
    if (/^You reached the hourly keyless Tavily limit|^Retry after:/im.test(content) ||
        (result.isError && /rate limit/i.test(content))) {
      console.log(content);
      break; // 遵守服务端限流，当前批次不继续发送请求。
    }
    // 全文只存原始文件，控制台保留题名、URL与短摘要，避免双栏论文挤占输出。
    if (mode === 'extract') console.log(content.slice(0,20000));
    else console.log(content.split(/(?=^Title: )/m).map(x=>x.split('Raw Content:')[0].slice(0,1800)).join('\n'));
    }
  }
}
