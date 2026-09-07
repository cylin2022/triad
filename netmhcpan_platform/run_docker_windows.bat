@echo off
chcp 65001 >nul
echo ============================================================
echo   TRIAD Platform - Windows Docker Launcher
echo ============================================================
echo.

docker info >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Docker Desktop 尚未啟動，請先開啟 Docker Desktop 後重試！
    pause
    exit /b 1
)

echo [*] 正在檢查舊容器狀態...
docker stop triad_app >nul 2>&1
docker rm triad_app >nul 2>&1
docker stop netmhcpan-app >nul 2>&1
docker rm netmhcpan-app >nul 2>&1

docker image inspect triad-platform:v1.1 >nul 2>&1
if errorlevel 1 (
    docker image inspect netmhcpan-platform:v1.1 >nul 2>&1
    if not errorlevel 1 (
        echo [*] 正在建立映像檔標籤 triad-platform:v1.1...
        docker tag netmhcpan-platform:v1.1 triad-platform:v1.1
    )
)

set MOUNT_ARG=
if exist "%~dp0netMHCpan-4.2" (
    echo [*] 偵測到本地 netMHCpan-4.2 資料夾，自動掛載...
    set MOUNT_ARG=-v "%~dp0netMHCpan-4.2":/opt/netMHCpan-4.2
) else if exist "%~dp0..\netMHCpan-4.2" (
    echo [*] 偵測到父層 netMHCpan-4.2 資料夾，自動掛載...
    set MOUNT_ARG=-v "%~dp0..\netMHCpan-4.2":/opt/netMHCpan-4.2
)

echo [*] 正在啟動 TRIAD 平台容器 (triad_app)...
docker run -d --name triad_app -p 5001:5001 %MOUNT_ARG% triad-platform:v1.1

if errorlevel 1 (
    echo.
    echo [ERROR] 啟動失敗！請確認映像檔已載入。
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

start http://localhost:5001

pause
