function problem_data = generateProblem(geometry, force, boundaries)

switch geometry
    case "square"
        problem_data = square_shell_problem(force, boundaries);
    case "hemi"
        problem_data = hemispherical_shell_problem(10,10);
end

end