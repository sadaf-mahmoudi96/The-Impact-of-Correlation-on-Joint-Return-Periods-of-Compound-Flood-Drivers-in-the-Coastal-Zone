To install NSVineCopula:

# 1. Go here:
https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html
# 2. Install Rtools42.
Restart RStudio.
Check if R can find make
Sys.which("make")
"C:\\rtools42\\usr\\bin\\make.exe"

install.packages(
  "https://cran.r-project.org/src/contrib/Archive/CDVineCopulaConditional/CDVineCopulaConditional_0.1.1.tar.gz",
  repos = NULL,
  type = "source"
)

# 3. install:
install.packages(
"https://github.com/agopb/NSVineCopula/blob/main/NSVineCopula_1.0.0.tar.gz",
repos = NULL,
type = "source"
)

OR

# 1) Download NSVineCopula source package from GitHub
url <- "https://github.com/agopb/NSVineCopula/raw/main/NSVineCopula_1.0.0.tar.gz"
tarfile <- file.path(tempdir(), "NSVineCopula_1.0.0.tar.gz")

download.file(url, tarfile, mode = "wb")

# 2) Unpack into a fresh folder
pkg_parent <- file.path(tempdir(), "NSVineCopula_patch_pointer")

if (dir.exists(pkg_parent)) {
  unlink(pkg_parent, recursive = TRUE)
}

dir.create(pkg_parent)

untar(tarfile, exdir = pkg_parent)

pkg_dir <- file.path(pkg_parent, "NSVineCopula")

# 3) Patch ALL .c files
c_files <- list.files(
  file.path(pkg_dir, "src"),
  pattern = "\\.c$",
  full.names = TRUE
)

for (cfile in c_files) {
  x <- readLines(cfile, warn = FALSE)
 
  # Add standard C library header for calloc/free
  if (!any(grepl("#include\\s*<stdlib\\.h>", x))) {
    x <- c("#include <stdlib.h>", x)
  }
 
  # Replace Calloc(n, double), Calloc(n, double*), Calloc(n, int), Calloc(n, int*)
  x <- gsub(
    "Calloc\\(([^,]+),\\s*double\\s*\\*\\s*\\)",
    "calloc(\\1, sizeof(double*))",
    x
  )
 
  x <- gsub(
    "Calloc\\(([^,]+),\\s*int\\s*\\*\\s*\\)",
    "calloc(\\1, sizeof(int*))",
    x
  )
 
  x <- gsub(
    "Calloc\\(([^,]+),\\s*double\\s*\\)",
    "calloc(\\1, sizeof(double))",
    x
  )
 
  x <- gsub(
    "Calloc\\(([^,]+),\\s*int\\s*\\)",
    "calloc(\\1, sizeof(int))",
    x
  )
 
  # Replace Free(x) with free(x)
  x <- gsub(
    "\\bFree\\(([^\\)]+)\\)",
    "free(\\1)",
    x
  )
 
  writeLines(x, cfile)
}

# 4) Install patched package
install.packages(pkg_dir, repos = NULL, type = "source")
