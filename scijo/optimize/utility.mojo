

struct RootResults[dtype: DType = DType.float64]():
    var root: Scalar[dtype]
    var iterations: Int
    var function_calls: Int
    var converged: Bool
    var flag: String
    var method: String

    fn __init__(out self,
        root: Scalar[dtype],
        iterations: Int,
        function_calls: Int,
        converged: Bool,
        flag: String,
        method: String,
    ):
        self.root = root
        self.iterations = iterations
        self.function_calls = function_calls
        self.converged = converged
        self.flag = flag
        self.method = method

    fn __str__(self) raises -> String:
        return String(
            "RootResults(root={}, iterations={}, function_calls={}, "
            "converged={}, flag='{}', method='{}')"
        ).format(
            self.root,
            self.iterations,
            self.function_calls,
            self.converged,
            self.flag,
            self.method,
        )

    fn write_to[W: Writer](self, mut writer: W):
        try:
            writer.write(
                String(
                    "Root Results\n"
                    "============\n"
                    "Root          : {}\n"
                    "Iterations    : {}\n"
                    "Function Calls: {}\n"
                    "Converged     : {}\n"
                    "Flag          : {}\n"
                    "Method        : {}\n"
                ).format(
                    self.root,
                    self.iterations,
                    self.function_calls,
                    self.converged,
                    self.flag,
                    self.method,
                )
            )
        except e:
            writer.write("Error displaying RootResults: " + String(e) + "\n")
