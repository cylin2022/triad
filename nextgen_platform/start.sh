#!/bin/bash
# start.sh — Launch TRIAD Platform
#
# Usage:
#   ./start.sh          # start on port 5001
#   ./start.sh 8080     # start on custom port

PORT=${1:-5001}
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "🧬 TRIAD: T-Cell Epitope Ranking, Immunogenicity, and Allele Distribution Framework"
echo "   Starting on http://localhost:${PORT}"
echo "   Press Ctrl+C to stop"
echo ""

# Kill any existing instance on the port
fuser -k ${PORT}/tcp 2>/dev/null

cd "$SCRIPT_DIR"
FLASK_PORT=${PORT} python3 -c "
import os
os.environ.setdefault('FLASK_PORT', '${PORT}')
from app import app
app.run(debug=False, host='0.0.0.0', port=${PORT}, threaded=True, use_reloader=False)
"
