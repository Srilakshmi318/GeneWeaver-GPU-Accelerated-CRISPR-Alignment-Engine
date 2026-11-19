
import time
from textual.app import App, ComposeResult
from textual.containers import Vertical
from textual.widgets import Header, Footer, Static, ProgressBar


class GeneWeaverApp(App):

    TITLE = "GENEWEAVER - CRISPR GPU ENGINE"
    SUB_TITLE = "Genomic Analysis Dashboard"

    CSS = """
    Screen {
        align: center middle;
    }

    #dashboard {
        width: 80%;
        height: auto;
        border: round green;
        padding: 2;
    }

    .info {
        color: cyan;
        margin: 1;
    }
    """

    def compose(self) -> ComposeResult:
        yield Header()

        with Vertical(id="dashboard"):
            yield Static(
                "GENEWEAVER - INTEGRATED CRISPR ENGINE",
                classes="info"
            )
            yield Static("GPU: Checking availability...", id="gpu")
            yield Static("Status: Ready", id="status")
            yield Static("CRISPR scanning dashboard", classes="info")
            yield ProgressBar(total=100, show_eta=False, id="progress")

        yield Footer()

    def on_mount(self):
        self.query_one("#gpu", Static).update(
            "GPU: Run CUDA scanning in Google Colab"
        )
        self.query_one("#status", Static).update(
            "Status: Dashboard started successfully"
        )


if __name__ == "__main__":
    GeneWeaverApp().run()
