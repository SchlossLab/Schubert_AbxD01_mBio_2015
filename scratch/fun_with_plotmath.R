my_names <- c("Name1", "Name2", "Name3", "Name4")
my_numbers <- c(2,4,5,9)
my_data <- c(24, 47, 55, 90)


par(mar=c(4,10,1,1))
plot(x=my_data, y=1:4, yaxt="n", ylab="", xlim=c(0,100))
)

# axis(2, label=paste0(my_names, " (Number ", my_numbers, ")"), at=1:4, las=1)

# no... just says "expresssion(...)"
# axis(2, label=paste0(expression(italic(my_names)), " (Number ", my_numbers, ")"), at=1:4, las=1)

# no... just italicizes "expresssion(...)"
#axis(2, label=parse(text=paste0(expression(italic(my_names)), "~(Number~", my_numbers, ")")), at=1:4, las=1)


# no... expression makes vector into a single value
#axis(2, label=expression(paste0(italic(my_names), " (Number ", my_numbers, ")")), at=1:4, las=1)

my_text <- paste0("italic(", my_names, ")~(Number~", my_numbers, ")")
axis(2, label=parse(text=my_text), at=1:4, las=1)


m <- "2.3~x~10^6"
l <- "1.2~x~10^6"
u <- "4.5~x~10^6"

my_string <- paste0(m, "~(", l, "~-~", u, ")")
text(x=50, y=2.5, parse(text=my_string))
