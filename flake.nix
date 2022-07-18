{
  inputs = {
    nixpkgs = {
      url = "github:nixos/nixpkgs/nixos-unstable";
    };

    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = import nixpkgs { inherit system; };
      in rec {
        packages.apm = pkgs.stdenv.mkDerivation {
          pname = "apm";
          version = "0.1.0";
          src = self;

          meta = with nixpkgs.lib; {
            description = "Anytime Performance Model";
            license = licenses.mit;
          };

          # Build dependencies
          nativeBuildInputs = with pkgs; [
            cmake
            ninja
            doxygen
            catch2
            gbenchmark
          ];
        };

        packages.default = self.packages.${system}.apm;
      });
}
