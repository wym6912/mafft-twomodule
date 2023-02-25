# MAFFT Module

This is the modified part of `mafft` modified by [gitee@wym6912](https://gitee.com/wym6912) / [github@wym6912](https://github.com/wym6912), which is a module of the other program. It has five subprograms: `staralign`, `profilealign`, `profile_two_align`, `fragalign` and `findsim`.

## How to compile

1. Compile with multi-threading support in Linux/WSL (default):  `make -j16`
2. Compile in Windows: 
   - required: `Visual Studio 2022`
   - Open `subwmsa.sln`, then choose `Build` -> `Build Solution`.
   - Open `x64\Debug` folder, you will find five files: `staralign.exe`, `profilealign.exe`, `fragalign.exe`, `findsim.exe` and `pthread.dll`. 
     - IMPORTANT: `pthread.dll` is necessary for other exe files.
     - Note: If you want to release this solution, please switch it to `Debug` mode. Select `Properties`, choose `Configuration Properties`, then press `Configuration Manager...` button, change `Active Solution configuration` to `Release`.

## How to use these programs

You can see `*-command.txt` for more infomation.

