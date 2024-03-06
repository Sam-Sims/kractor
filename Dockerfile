#https://dev.to/rogertorres/first-steps-with-docker-rust-30oi
FROM rust:1.76.0-slim as build
LABEL authors="sammy"

# Create the empty cargo project
RUN USER=root cargo new --bin kractor
workdir /kractor

# Copy the Cargo.toml and Cargo.lock files
COPY ./Cargo.toml ./Cargo.toml
COPY ./Cargo.lock ./Cargo.lock

# Build the dependencies
RUN cargo build --release
RUN rm src/*.rs
COPY ./src ./src

# Build
RUN rm ./target/release/deps/kractor*
RUN cargo build --release

FROM rust:1.76.0-slim
COPY --from=build /kractor/target/release/kractor /usr/local/bin/kractor

ENTRYPOINT ["./usr/local/bin/kractor"]