#!/bin/sh
# entrypoint.sh — checks netMHCpan-4.2 install status before launching the platform.
# Missing binary/data is NOT a fatal error: app.py already serves a /setup
# wizard so the user can install it from the browser after the container
# is up.
set -e

BIN="${NETMHC_HOME:-/opt/netMHCpan-4.2}/Linux_x86_64/bin/netMHCpan-4.2"

echo "============================================================"
echo " TRIAD Platform (port ${FLASK_PORT:-5001}) — v1.1"
echo "============================================================"
if [ -x "$BIN" ]; then
    echo "[entrypoint] netMHCpan-4.2 found at: $BIN"
else
    echo "[entrypoint] netMHCpan-4.2 NOT found at: $BIN"
    echo "[entrypoint] Starting in setup mode — open http://localhost:${FLASK_PORT:-5001}/setup"
    echo "[entrypoint] to install it (mount a volume at \$NETMHC_HOME, or upload the"
    echo "[entrypoint] tar.gz you downloaded from DTU through the web wizard)."
fi
echo "============================================================"

exec "$@"
