import http.server
import socketserver
from threading import Timer
import webbrowser


def is_port_in_use(port):
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

def check_port():
    for i in range(8000, 8100):
            if is_port_in_use(i) == False:
                port = i
                break
    return i

def run_server(PORT):
    Handler = http.server.SimpleHTTPRequestHandler
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        print("serving at port", PORT)
        httpd.serve_forever()

def open_localhost(port, path):
    webbrowser.open_new(f"http://localhost:{port}/{path}")

def run_sunburst(path='sunburst/'):
    port = check_port()
    Timer(1, open_localhost, args=[port, path]).start();
    run_server(port)
