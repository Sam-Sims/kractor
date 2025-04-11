use chrono::Local;
use clap::Parser;
use color_eyre::eyre::bail;
use color_eyre::Result;
use env_logger::fmt::Color;
use env_logger::Builder;
use log::LevelFilter;
use std::io::Write;

pub mod extract;
pub mod parsers;
pub use crate::cli::Cli;
pub mod cli;
pub mod kractor;

use kractor::Kractor;

fn init_logging(verbose: bool) {
    let level_filter = if verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };

    Builder::new()
        .format(|buf, record| {
            let mut style = buf.style();
            style.set_color(match record.level() {
                log::Level::Trace => Color::Magenta,
                log::Level::Debug => Color::Blue,
                log::Level::Info => Color::Green,
                log::Level::Warn => Color::Yellow,
                log::Level::Error => Color::Red,
            });

            writeln!(
                buf,
                "{} [{}] - {}",
                Local::now().format("[%H:%M:%S]"),
                style.value(record.level()),
                record.args()
            )
        })
        .filter(None, level_filter)
        .init();
}

fn main() -> Result<()> {
    let args = Cli::parse();
    init_logging(args.verbose);

    color_eyre::install()?;

    if args.input.len() != args.output.len() {
        bail!("Number of input and output files must match");
    }

    let mut app = Kractor::new(args);
    app.run()?;

    Ok(())
}
