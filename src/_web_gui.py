import json
import os
import uuid

# jscircuit: a browser-based circuit editor producing netlists in the
# same format as qucat's native tkinter GUI (see qucat.GUI).
# https://github.com/qucat/jscircuit
JSCIRCUIT_URL = "https://qucat.github.io/jscircuit/app/jscircuit.html"


def web_gui(filename=None, width=700, height=500):
    r'''Opens a browser-based circuit editor inside a Jupyter notebook cell.

    This is an alternative to :class:`qucat.GUI` for situations where the
    native tkinter interface is unavailable or impractical, for example
    when working on a remote server or inside JupyterHub. It embeds
    `jscircuit <https://github.com/qucat/jscircuit>`_, a browser-based
    circuit editor, together with a small text editor used to save the
    resulting netlist to a file.

    Parameters
    ----------
    filename:   string, optional
                Path to the file the netlist will be saved to. If the
                file already exists, its contents are loaded into the
                text editor, and the circuit is also loaded into the
                jscircuit editor itself. If not given, the file name
                can instead be typed into the text editor once the
                cell has been run.
    width:      int, optional
                Width in pixels of the jscircuit editor (default 700).
                The editor can also be resized by dragging its
                bottom-right corner.
    height:     int, optional
                Height in pixels of the jscircuit editor (default 500).

    Returns
    -------
    None

    Notes
    -----

    This function does not return a :class:`qucat.Qcircuit`: since the
    circuit is drawn in the browser, the netlist only becomes available
    once the user has saved it, which happens after this function
    returns. The typical workflow is thus split over two cells:

    1. Call ``qucat.web_gui(filename)`` in one cell, draw the circuit in
       the jscircuit editor, copy the resulting netlist, paste it into
       the text box below the editor, then click "Save".
    2. In a later cell, load the saved netlist with
       ``circuit = qucat.GUI(filename, edit=False)``.

    Loading an existing netlist into the jscircuit editor relies on
    jscircuit only accepting such requests from a same-origin,
    ``localhost``, or ``*.github.io`` parent page (a restriction of
    jscircuit itself, for security). This holds for a notebook running
    on your own machine, but on a remote JupyterHub reached through a
    different domain, the circuit will not auto-load: the netlist is
    still pre-filled in the text box below the editor, from where it
    can be pasted into jscircuit manually (Ctrl+V).

    Requires the optional ``ipywidgets`` dependency, install it with
    ``pip install ipywidgets`` or ``pip install qucat[notebook]``.
    '''
    try:
        import ipywidgets as widgets
        from IPython.display import display, HTML
    except ImportError as e:
        raise ImportError(
            "qucat.web_gui requires the 'ipywidgets' package, which is "
            "not installed.\nInstall it with:\n\n"
            "    pip install ipywidgets\n\n"
            "or install qucat with the 'notebook' extra:\n\n"
            "    pip install qucat[notebook]"
        ) from e

    # Disable the native drag-to-resize handle on the textarea, it is
    # redundant with the resize handle of the jscircuit editor below.
    display(HTML("<style>textarea { resize: none !important; }</style>"))

    initial_text = ""
    if filename is not None and os.path.isfile(filename):
        with open(filename, "r", encoding="utf-8") as f:
            initial_text = f.read()

    container_id = "qucat-web-gui-%s" % uuid.uuid4().hex

    # jscircuit posts {type: 'appReady'} to its parent window once loaded,
    # and after that accepts {type: 'loadCircuit', netlist: <string>} to
    # populate the editor. We instead trigger on the iframe's own 'load'
    # event, which (since main.js is a module script) only fires after
    # jscircuit has already sent 'appReady', so the app is guaranteed to
    # be ready to receive the netlist by then.
    autoload_script = ""
    if initial_text:
        # Escape "</" so the embedded JSON can't prematurely close the
        # surrounding <script> tag if the netlist ever contained it.
        netlist_json = json.dumps(initial_text).replace("</", "<\\/")
        autoload_script = '''
        <script>
        (function() {
            var iframe = document.getElementById("%s").querySelector("iframe");
            iframe.addEventListener("load", function() {
                iframe.contentWindow.postMessage(
                    {type: "loadCircuit", netlist: %s}, "*");
            });
        })();
        </script>
        ''' % (container_id, netlist_json)

    display(HTML('''
    <div id="{container_id}" style="
        resize: both;
        overflow: hidden;
        width: {width}px;
        height: {height}px;
        min-width: 200px;
        min-height: 150px;
        border: 1px solid #ccc;
        box-sizing: border-box;
    ">
      <iframe
        src="{url}"
        allow="clipboard-write; clipboard-read"
        width="100%" height="100%"
        style="border: none; display: block;">
      </iframe>
    </div>
    {autoload_script}
    '''.format(container_id=container_id, width=width, height=height,
               url=JSCIRCUIT_URL, autoload_script=autoload_script)))

    text_area = widgets.Textarea(
        value=initial_text,
        placeholder="Paste the netlist copied from jscircuit here...",
        layout=widgets.Layout(width="500px", height="200px"),
    )

    filename_input = widgets.Text(
        value=filename or "",
        placeholder="e.g. netlist.txt",
        description="File name:",
        layout=widgets.Layout(width="300px"),
    )

    save_button = widgets.Button(
        description="Save",
        button_style="primary",
        icon="save",
    )

    status_label = widgets.Label(value="")

    def on_save_clicked(b):
        fname = filename_input.value.strip()
        content = text_area.value

        if not fname:
            status_label.value = "Please enter a file name."
            return

        try:
            directory = os.path.dirname(fname)
            if directory:
                os.makedirs(directory, exist_ok=True)
            with open(fname, "w", encoding="utf-8") as f:
                f.write(content)
            status_label.value = "Saved to '%s' (%d chars)" % (fname, len(content))
        except Exception as e:
            status_label.value = "Error: %s" % e

    save_button.on_click(on_save_clicked)

    ui = widgets.VBox([
        widgets.Label("Netlist"),
        text_area,
        widgets.HBox([filename_input, save_button]),
        status_label,
    ])

    display(ui)
