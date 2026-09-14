import ast


def load_functions(script_path, extra_namespace=None):
    """
    Load all top-level functions from a Python script
    without executing top-level code such as sys.argv main calls.
    """

    with open(script_path, "r") as f:
        source = f.read()

    tree = ast.parse(source)

    function_nodes = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
    ]

    module = ast.Module(
        body=function_nodes,
        type_ignores=[]
    )

    namespace = {}

    if extra_namespace is not None:
        namespace.update(extra_namespace)

    exec(
        compile(module, script_path, "exec"),
        namespace
    )

    return namespace