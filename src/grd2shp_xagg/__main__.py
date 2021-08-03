"""Command-line interface."""
import click


@click.command()
@click.version_option()
def main() -> None:
    """Grd2Shp_Xagg."""


if __name__ == "__main__":
    main(prog_name="grd2shp_xagg")  # pragma: no cover
