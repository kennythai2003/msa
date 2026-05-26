from cffi import FFI
import os
import subprocess

# help from Claude for cffi set up
ffi = FFI()

ffi.cdef("""
    int msa_align3(const char *u, int m,
                   const char *v, int n,
                   const char *w, int k,
                   int *aln_len, int *perfect,
                   char *out_u, char *out_v, char *out_w);
""")

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "align3_core.cpp")) as f:
    source = f.read()

try:
    sdk = subprocess.check_output(
        ["xcrun", "--show-sdk-path"], stderr=subprocess.DEVNULL
    ).decode().strip()
    extra_args = ["-O3", f"-isysroot{sdk}"]
except Exception:
    extra_args = ["-O3"]

ffi.set_source("_align3_core", source,
               source_extension=".cpp",
               extra_compile_args=extra_args,
               extra_link_args=extra_args)

if __name__ == "__main__":
    ffi.compile(tmpdir=here, verbose=True)
