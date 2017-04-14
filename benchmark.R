

true_beta1 = seq(0, 2, 0.4)
true_beta0 = seq(0, 2, 0.4)

nb_sample = seq(50, 1000, 150)

nb_iter = seq(100, 1000,200)



iter = 0

for (i in true_beta0) {
  for(j in true_beta1) {
    for(k in nb_sample) {
      for(t in nb_iter) {
        iter = iter +1
      }
    }
  }
}
print(iter)