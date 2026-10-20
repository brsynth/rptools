DEFAULT_DELIMITER = ","


def add_arguments(parser):
    parser.add_argument(
        "--pathways",
        required=True,
        type=str,
        nargs="+",
        help="Pathways (rpSBML) to rank",
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        default=DEFAULT_DELIMITER,
        help=f"Set delimiter character for output (default: {DEFAULT_DELIMITER})",
    )
    return parser
