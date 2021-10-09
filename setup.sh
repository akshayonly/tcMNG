mkdir -p ~/.streamlit/

echo "[theme]
primaryColor='#5bd0f9'
backgroundColor='#222222'
secondaryBackgroundColor='#4b4b4b'
textColor='#f1f2f6'
[server]
headless = true
port = $PORT
enableCORS = false
" > ~/.streamlit/config.toml