### Troubleshooting

#### Installing popt

Most model Linux systems don't include the popt commandline parsing C library anymore. In case you should encounter a compilation error referencing a missing popt library, you can follow the steps below to install popt:

1. Download popt and un-tar the compressed archive file:

  wget http://rpm5.org/files/popt/popt-1.16.tar.gz popt-1.16.tar.gz
  tar -xvzf popt-1.16.tar.gz
  
2. Configure the popt installer:

  cd popt-1.16
  ./configure --prefix=/usr --disable-static &&
  make
  
3. Install popt:

- Note that this step may need to be performed as root.

  make install