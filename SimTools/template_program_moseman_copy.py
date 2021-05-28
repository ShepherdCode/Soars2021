import os, sys, traceback, argparse

class Template_Class():
    '''How to write a class.'''
    def __init__(self,debug=False):
        '''How to write a constructor.'''
        self.debug=debug
    def show_off(self):
        '''How to access instance variables.'''
        print("TemplateClass")
        if self.debug:
            print("\tIn debug mode.")
    def write_file(self,filename,lines=10):
        '''How to set parameters with defaults.'''
        if self.debug:
            print("\tWriting %d lines to file: %s."%
                (lines,filename))
        with open(filename, 'w') as outfa:
            for line in range(0,lines):
                outfa.write("2 ^ " + str(line) + " = " + str(2**line) + "\n")
        print("\tDon't forget to delete the file!")
    def no_op(self):
        '''How to write a method that does nothing.'''
        pass

def args_parse():
    '''How to parse command-line arguments.'''
    global args
    parser = argparse.ArgumentParser(
        description='Bare bones Python program.')
    parser.add_argument(
        'numlines',
        help='output file size (10)',
        type=int)
    parser.add_argument(
        'outfile',
        help='output filename (fasta)',
        type=str)
    parser.add_argument(
        '--debug',
        help='Print traceback after exception.',
        action='store_true')
    args = parser.parse_args()

if __name__ == "__main__":
    '''How to start a program from the command line.'''
    try:
        args_parse()
        numlines=args.numlines
        outfile=args.outfile
        debug=args.debug
        tmp = Template_Class(debug)
        tmp.show_off()
        tmp.write_file(outfile,numlines)
        tmp.no_op()
    except Exception:
        print()
        if args.debug:
            print(traceback.format_exc())
        else:
            print("There was an error.")
            print("Run with --debug for traceback.")
