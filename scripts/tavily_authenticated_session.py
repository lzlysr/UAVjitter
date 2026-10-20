"""以隐藏输入建立临时 Tavily 会话；密钥不落盘、不放入命令行参数。"""
import getpass
import os
import subprocess
from pathlib import Path

if __name__ == "__main__":
    session_env = os.environ.copy()
    session_env["TAVILY_API_KEY"] = getpass.getpass("请输入 Tavily API key（隐藏输入，仅本次会话使用）: ")
    if not session_env["TAVILY_API_KEY"]:
        raise SystemExit("未提供密钥。")
    script = Path(__file__).with_name("tavily_literature_search.mjs")
    try:
        subprocess.run(["node", str(script), "session"], env=session_env, check=True)
    except KeyboardInterrupt:
        print("已结束临时认证会话。")
    finally:
        session_env.pop("TAVILY_API_KEY", None)
