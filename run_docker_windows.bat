@echo off
chcp 65001 >nul
echo ============================================================
echo   TRIAD Platform - Windows Docker Launcher
echo ============================================================
echo.

:: 1. 檢查 Docker 是否運行
docker info >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Docker Desktop 尚未啟動，請先開啟 Docker Desktop 後重試！
    pause
    exit /b 1
)

:: 2. 清理舊容器 (triad_app 或舊名稱)
echo [*] 正在檢查舊容器狀態...
docker stop triad_app >nul 2>&1
docker rm triad_app >nul 2>&1
docker stop netmhcpan-app >nul 2>&1
docker rm netmhcpan-app >nul 2>&1

:: 3. 檢查映像檔 (若有舊版名稱則自動別名標記為 triad-platform:v1.1)
docker image inspect triad-platform:v1.1 >nul 2>&1
if errorlevel 1 (
    docker image inspect netmhcpan-platform:v1.1 >nul 2>&1
    if not errorlevel 1 (
        echo [*] 發現映像檔 netmhcpan-platform:v1.1，正在建立標籤 triad-platform:v1.1...
        docker tag netmhcpan-platform:v1.1 triad-platform:v1.1
    ) else (
        if exist "%~dp0netmhcpan-platform-v1.1.tar" (
            echo [*] 正在載入映像檔 netmhcpan-platform-v1.1.tar...
            docker load -i "%~dp0netmhcpan-platform-v1.1.tar"
            docker tag netmhcpan-platform:v1.1 triad-platform:v1.1 >nul 2>&1
        ) else (
            echo [ERROR] 找不到映像檔 triad-platform:v1.1！
            echo 請確認 Docker Desktop 已載入映像檔。
            pause
            exit /b 1
        )
    )
)

:: 4. 檢查本地 netMHCpan-4.2 資料夾以掛載 Volume
set MOUNT_ARG=
if exist "%~dp0netMHCpan-4.2" (
    echo [*] 偵測到本地 netMHCpan-4.2 資料夾，自動掛載...
    set MOUNT_ARG=-v "%~dp0netMHCpan-4.2":/opt/netMHCpan-4.2
)

:: 5. 啟動 TRIAD 容器 (名稱: triad_app，連接埠: 5001)
echo [*] 正在啟動 TRIAD 平台容器 (triad_app)...
docker run -d --name triad_app -p 5001:5001 %MOUNT_ARG% triad-platform:v1.1

if errorlevel 1 (
    echo.
    echo [ERROR] 啟動失敗！
    pause
    exit /b 1
)

echo.
echo ============================================================
echo   TRIAD 平台啟動成功！
echo   容器名稱: triad_app
echo   請在瀏覽器開啟: http://localhost:5001
echo   (注意: Windows 請勿使用 0.0.0.0:5001)
echo ============================================================
echo.

:: 6. 自動喚起預設瀏覽器開啟 TRIAD 網頁
start http://localhost:5001

pause
