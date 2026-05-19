{
  pkgs,
  lib,
  config,
  inputs,
  ...
}:
{
  languages.python.uv = {
    enable = true;
    sync.enable = true;
  };
  env.LD_LIBRARY_PATH = "${lib.makeLibraryPath [
    pkgs.stdenv.cc.cc.lib
    pkgs.zlib
  ]}";
}
