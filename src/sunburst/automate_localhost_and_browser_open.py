'''
import tkinter as tk
from tkinterhtml import HtmlFrame
from tk_html_widgets import HTMLLabel


root = tk.Tk()
root.title("Sunburst plotting")
screen_width = root.winfo_screenwidth()*0.75
screen_height = root.winfo_screenheight()*0.75
print(screen_width)
print(type(screen_width))
root.geometry(f'{int(screen_width)}x{int(screen_height)}')#widthxheight


html_label = HTMLLabel(root, html='<h1 style="color: red; text-align: center"> Hello World </H1>')
html_label.pack(fill="both", expand=True)
html_label.fit_height()
'''
from threading import Timer
import webbrowser
def run_server():
    try:
        # Python 2
        from SimpleHTTPServer import test, SimpleHTTPRequestHandler
    except ImportError:
        # Python 3
        from http.server import test, SimpleHTTPRequestHandler
        test(SimpleHTTPRequestHandler)


def open_localhost():
    webbrowser.open_new("http://localhost:8000/sunburst/adjusted_sunburst.html")


def run_sunburst():
    Timer(1, open_localhost).start();
    run_server()

run_sunburst()
