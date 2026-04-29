%%
clc;
clear;

%%
load("results/cfb_coeffs.mat")

precision = 17;
output_path = fullfile("results", "raipal_solution", "cfbSolution.hpp");
output_dir = fileparts(output_path);

if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

fid = fopen(output_path, "w");
if fid < 0
    error("Could not open %s for writing.", output_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, "#pragma once\n\n");
fprintf(fid, "#include <array>\n");
fprintf(fid, "#include <cstddef>\n\n");
fprintf(fid, "// Generated from results/cfb_coeffs.mat\n");
fprintf(fid, "// Coefficients are highest to lowest degree for normalized x.\n");
fprintf(fid, "// x = 2.0 * (q - min) / (max - min) - 1.0;\n\n");
fprintf(fid, "namespace cfb {\n\n");

fprintf(fid, "template <std::size_t N>\n");
fprintf(fid, "struct cfbParam {\n");
fprintf(fid, "    double min;\n");
fprintf(fid, "    double max;\n");
fprintf(fid, "    std::array<double, N> coeff0;\n");
fprintf(fid, "    std::array<double, N> coeff1;\n");
fprintf(fid, "    std::array<double, N> coeff2;\n");
fprintf(fid, "};\n\n");

if exist("coeffs_b", "var")
    print_cfb_param(fid, "fromJoint", coeffs_b.to_a, coeffs_b.to_da, coeffs_b.to_dda, precision);
end

if exist("coeffs_a", "var")
    print_cfb_param(fid, "fromActuator", coeffs_a.to_b, coeffs_a.to_db, coeffs_a.to_ddb, precision);
end

fprintf(fid, "inline constexpr double actuatorInertia = 0.00134589; // corresponds to the inertia of elb_input / elb_bar\n\n");

fprintf(fid, "template <std::size_t N>\n");
fprintf(fid, "constexpr double evalPoly(const std::array<double, N>& p, double x) {\n");
fprintf(fid, "    double y = p[0];\n");
fprintf(fid, "    for (std::size_t i = 1; i < N; ++i) {\n");
fprintf(fid, "        y = y * x + p[i];\n");
fprintf(fid, "    }\n");
fprintf(fid, "    return y;\n");
fprintf(fid, "}\n\n");

fprintf(fid, "template <std::size_t N>\n");
fprintf(fid, "constexpr void evalCfb(const cfbParam<N>& param, double q, double& value, double& first, double& second) {\n");
fprintf(fid, "    const double x = 2.0 * (q - param.min) / (param.max - param.min) - 1.0;\n");
fprintf(fid, "    value = evalPoly(param.coeff0, x);\n");
fprintf(fid, "    first = evalPoly(param.coeff1, x);\n");
fprintf(fid, "    second = evalPoly(param.coeff2, x);\n");
fprintf(fid, "}\n\n");

fprintf(fid, "} // namespace cfb\n");

fprintf("Saved C++ header to %s\n", output_path);

function print_cfb_param(fid, var_name, coeff0, coeff1, coeff2, precision)
    validate_compatible_coeffs(var_name, coeff0, coeff1, coeff2);
    n = numel(coeff0.coeffs);

    fprintf(fid, "inline constexpr cfbParam<%d> %s = {\n", n, var_name);
    fprintf(fid, "    %s,\n", format_scalar(coeff0.min, precision));
    fprintf(fid, "    %s,\n", format_scalar(coeff0.max, precision));
    print_vector_expr(fid, "coeff0", coeff0.coeffs(:), precision, true);
    print_vector_expr(fid, "coeff1", coeff1.coeffs(:), precision, true);
    print_vector_expr(fid, "coeff2", coeff2.coeffs(:), precision, false);
    fprintf(fid, "};\n\n");
end

function validate_compatible_coeffs(var_name, coeff0, coeff1, coeff2)
    if coeff0.degree ~= coeff1.degree || coeff0.degree ~= coeff2.degree
        error("%s has mismatched polynomial degrees.", var_name);
    end

    if coeff0.min ~= coeff1.min || coeff0.min ~= coeff2.min || ...
       coeff0.max ~= coeff1.max || coeff0.max ~= coeff2.max
        error("%s has mismatched normalization ranges.", var_name);
    end
end

function print_vector_expr(fid, label, coeffs, precision, trailing_comma)
    n = numel(coeffs);

    fprintf(fid, "    std::array<double, %d>{\n", n);
    for i = 1:4:n
        line_end = min(i + 3, n);
        fprintf(fid, "        ");

        for j = i:line_end
            fprintf(fid, "%s", format_scalar(coeffs(j), precision));

            if j < n
                fprintf(fid, ", ");
            end
        end

        fprintf(fid, "\n");
    end

    if trailing_comma
        fprintf(fid, "    }, // %s\n", label);
    else
        fprintf(fid, "    }  // %s\n", label);
    end
end

function str = format_scalar(value, precision)
    str = sprintf("%.*g", precision, value);

    if ~contains(str, ".") && ~contains(str, "e") && ~contains(str, "E")
        str = str + ".0";
    end
end
