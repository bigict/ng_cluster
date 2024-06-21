from ArgumentParserBuilder import NGArgumentParserBuilder


class ClusterArgumentParser(NGArgumentParserBuilder):
    def __init__(self):
        super().__init__()

        # Program description
        self.description = 'This is an argument parser for Cluster Analysis.'
        self.name = 'Cluster Analysis'
        self.version = '1.0'
        