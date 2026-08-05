version 1.0

workflow miniwdl_smoke {
  call write_message

  output {
    File result = write_message.result
  }
}

task write_message {
  command <<<
    printf 'shared miniwdl works\n' > result.txt
  >>>

  output {
    File result = "result.txt"
  }

  runtime {
    cpu: 1
    memory: "512 MiB"
    docker: "ubuntu:20.04"
  }
}
