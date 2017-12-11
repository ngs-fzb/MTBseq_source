/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in -- template for config.h.  If you use configure, this file
   provides #defines reflecting your configuration choices.  If you don't
   run configure, suitable conservative defaults will be used.

   This template file can be updated with autoheader, but do so carefully
   as autoheader adds #defines such as PACKAGE_* that we don't want.  */

/* Define to 1 if a SysV or X/Open compatible Curses library is present */
#define HAVE_CURSES 1

/* Define to 1 if library supports color (enhanced functions) */
#define HAVE_CURSES_COLOR 1

/* Define to 1 if library supports X/Open Enhanced functions */
/* #undef HAVE_CURSES_ENHANCED */

/* Define to 1 if <curses.h> is present */
/* #undef HAVE_CURSES_H */

/* Define to 1 if library supports certain obsolete features */
#define HAVE_CURSES_OBSOLETE 1

/* Define to 1 if the Ncurses library is present */
#define HAVE_NCURSES 1

/* Define to 1 if the NcursesW library is present */
/* #undef HAVE_NCURSESW */

/* Define to 1 if <ncursesw/curses.h> is present */
/* #undef HAVE_NCURSESW_CURSES_H */

/* Define to 1 if <ncursesw.h> is present */
/* #undef HAVE_NCURSESW_H */

/* Define to 1 if <ncurses/curses.h> is present */
/* #undef HAVE_NCURSES_CURSES_H */

/* Define to 1 if <ncurses.h> is present */
#define HAVE_NCURSES_H 1

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */
