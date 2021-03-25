rule test:
    input: "A.txt"

    resources:
        cpus = 16,
        time = "10:00",
        mem = "16G"

    params:
        partition = "talon"

    script: "test.R"