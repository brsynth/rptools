def add_arguments(parser):
    parser.add_argument(
        "--pathways",
        required=True,
        type=str,
        nargs="+",
        help="Pathways (rpSBML) to make statistics on",
    )
    return parser
