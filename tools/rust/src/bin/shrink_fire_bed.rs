fn main() {
    if let Err(error) = longpipe_tools::shrink_fire_bed::run() {
        eprintln!("shrink_fire_bed: {error:#}");
        std::process::exit(1);
    }
}
