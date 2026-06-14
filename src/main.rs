pub mod cli;
pub mod extract;
pub mod kractor;
pub mod parsers;

use std::io::Write;

use chrono::Local;
use clap::Parser;
use color_eyre::{Result, eyre::bail};
use env_logger::{Builder, fmt::Color};
use log::LevelFilter;

pub use crate::cli::Cli;

fn main() -> Result<()> {
    color_eyre::install()?;

    let args = Cli::parse();
    init_logging(args.verbose);

    if args.input.len() != args.output.len() {
        bail!("Number of input and output files must match");
    }

    kractor::run(args)?;

    Ok(())
}

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
