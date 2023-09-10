using Gradus



m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 10_000.0, deg2rad(60), 0.0)
d = ThinDisc(Gradus.isco(m), 500.0)


rendergeodesics(m, x, d, 2000.0)
