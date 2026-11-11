fn main() {
    if let Err(error) = longpipe_tools::copy_pbwdl::run() {
        eprintln!("copy_pbWDL: {error:#}");
        std::process::exit(1);
    }
}
