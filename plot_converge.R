d = read.table('Uso_filter10_fast0_seed100_iter10000.log.txt', sep=' ' )


pdf('fast0_converge.pdf')
plot(d[,1], d[,2], pch=16)
lines(d[,1], d[,2])
dev.off()

pdf('fast0_converge.log.pdf')
plot(d[,1], d[,2], log='x', pch=16)
lines(d[,1], d[,2])
dev.off()

