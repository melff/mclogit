# Deparse into a single string
deparse0 <- function(formula) paste(trimws(deparse(formula)),collapse=" ")

# Concatenate two formulae
c_formulae <- function(formula,extra){
    formula.deparsed <- deparse0(formula)
    extra.deparsed <- sub("~","+",deparse0(extra)) 
    as.formula(paste(formula.deparsed,
                     extra.deparsed),
               env=environment(formula))
}

# Check if formula
is_formula <- function(x)inherits(x,"formula")

# Subtitute "|" with "+"
random2formula <- function(r) {
    formula.deparsed <- deparse0(r$formula)
    gf <- paste(r$groups,collapse="+")
    as.formula(paste(formula.deparsed,
                     gf,sep="+"),
               env=environment(r$formula))
}
