#!/bin/bash
# sqfoam installer — sets up the `sqfoam` command and writes ~/.sqfoam.conf.
# Assumes OpenFOAM-6 + LIGGGHTS + CFDEMcoupling + libsqCfdem are already built
# (see docs/INSTALL_STACK.md to build the stack itself).
set -e
DEST="${1:-$HOME/.local/sqfoam}"
mkdir -p "$DEST"
cp -r sqfoam "$DEST/"
BIN="$HOME/.local/bin"; mkdir -p "$BIN"
cat > "$BIN/sqfoam" << SH
#!/bin/bash
exec python3 "$DEST/sqfoam/__main__.py" "\$@"
SH
chmod +x "$BIN/sqfoam"
echo "installed sqfoam -> $BIN/sqfoam"
case ":$PATH:" in *":$BIN:"*) : ;; *)
  echo "NOTE: add ~/.local/bin to PATH:  echo 'export PATH=\$HOME/.local/bin:\$PATH' >> ~/.bashrc" ;;
esac
echo "next:  sqfoam doctor"
