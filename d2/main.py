"""Main CLI application for the d2 network system.

Following the "One File Per Command" pattern from Typer documentation.
"""

import typer

# Import command modules
from d2.commands import (
    add,
    calc,
    collect,
    deprotonate,
    process_md,
    protonate,
    run_md,
    setup,
    tautomerize,
    ts_search,
    vis,
)

app = typer.Typer()

# Add all command modules to the main app
app.add_typer(setup.app, name=None)  # Top level commands
app.add_typer(vis.app, name=None)
app.add_typer(collect.app, name=None)
app.add_typer(calc.app, name=None)
app.add_typer(tautomerize.app, name=None)
app.add_typer(protonate.app, name=None)
app.add_typer(deprotonate.app, name=None)
# app.add_typer(show.app, name=None)
app.add_typer(add.app, name=None)
app.add_typer(run_md.app, name=None)
# app.add_typer(plot_reactions.app, name=None)
app.add_typer(process_md.app, name=None)
app.add_typer(ts_search.app, name=None)


if __name__ == "__main__":
    app()
