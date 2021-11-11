{
  description = "apm";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = import nixpkgs { inherit system; };
      in rec {
        packages.apm = pkgs.stdenv.mkDerivation {
          pname = "apm";
          version = "0.1.0";
          src = self;
          # Build dependencies
          nativeBuildInputs = with pkgs; [ meson ninja pkg-config ];
          # Run-time dependencies
          buildInputs = with pkgs; [ glpk gbenchmark ];
        };
        defaultPackage = self.packages.${system}.apm;
      });
}
