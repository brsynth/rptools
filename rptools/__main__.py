from os import path as os_path

def _cli():

    with open(os_path.join(os_path.dirname(os_path.abspath(__file__)), '.env'), 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('MODULES='):
                modules = line.splitlines()[0].split('=')[1].lower().split(',')

    print()
    print('Welcome to rpTools!')
    print()
    print('\'rptools\' is a package and cannot be directly executed. Executable tools are:')
    for module in modules:
        print('   - '+module)
    print()
    print('To find help for a specific tool, please type:')
    print('   python -m rptools.<tool_name> --help')
    print()

    return 0

import os
if __name__ == '__main__':
    _cli()
