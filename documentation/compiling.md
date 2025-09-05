# Compile Notes

Currently this project uses portable-simd which is a nightly only
feature. So, you have to use the nightly build to compile.

(rustup override set nightly)
(rustup override unset nightly)

Hopefully I'll be able to add in non-simd variants with optional
compile for either nightly or standard build. Or maybe the
portable-simd code will go mainstream soon. There are some other
compile notes in the top of the main.rs file. This also cross compiles
from Linux to windows or natively on windows.

If you look at the code, you can see the beginnings of my efforts to
get the non-simd stuff to work as an optional compile, but I'm
uploading in the current state because it works, and I want something
that works out in public.

