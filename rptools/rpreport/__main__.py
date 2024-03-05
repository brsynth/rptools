from rptools.rpreport.Args import add_arguments
from rptools import build_args_parser

from rptools.rpreport.rpreport import run_report


def entry_point():
    parser = build_args_parser(
        prog='rpreport',
        description='generates HTML pages to explore the main characteristics (thermodynamics, fluxes, number of '
                    'metabolic steps, reaction rule score) of pathways predicted with RetroPath suite',
        m_add_args=add_arguments
    )
    args = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    run_report(
        input_dir=args.input_dir,
        source_path=args.source_path,
        output_folder=args.output_folder,
        dev=args.dev,
        verbose=args.verbose,
        standalone=args.standalone,
    )


if __name__ == '__main__':
    entry_point()
