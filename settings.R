
# Increase terminal width in interactive sessions
terminal.width <- Sys.getenv('COLUMNS')
terminal.width <- ifelse(terminal.width == '', 80, terminal.width)
options(width=terminal.width)

# Read feature descriptions and easy_names
feature.desc.path <- file.path(project.dir, 'data-integration', 'feature-descriptions.txt')
desc.df <- read.delim(feature.desc.path)
feature.converter <- desc.df$easy_name
names(feature.converter) <- desc.df$feature

# Read pathophysiology
pathophys.path <- file.path(project.dir, 'data-integration', 'pathophysiology.txt')
pathophys.df <- read.delim(pathophys.path)

