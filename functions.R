
# Increase terminal width in interactive sessions
terminal.width <- Sys.getenv('COLUMNS')
terminal.width <- ifelse(terminal.width == '', 80, terminal.width)
options(width=terminal.width)



