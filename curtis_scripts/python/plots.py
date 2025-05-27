# Python plotting helper for CI

def save_plot(fig, name):
    """Save a matplotlib figure to the plots directory with a standard name."""
    import os
    outdir = os.path.join(os.path.dirname(__file__), '../../plots')
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=150)
