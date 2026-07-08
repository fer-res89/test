"""Small internal library shared by the forge and aeon CLIs.

Only truly generic, dependency-light helpers belong here (path/file
utilities). Tree-parsing helpers stay inside each tool's own package since
forge uses Bio.Phylo and aeon uses dendropy, and the two are not
currently interchangeable.
"""
