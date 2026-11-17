import numpy as np

def test_periodic_mapping():
    L = 1.0
    test_cases = [
        (-0.2, 0.8),
        (1.3, 0.3),
        (0.5, 0.5),
        (-1.1, 0.9),
        (2.7, 0.7),
    ]
    
    for x, expected in test_cases:
        result = (x + L) % L
        assert np.isclose(result, expected), f"Failed for {x}: {result} vs {expected}"
    print("所有测试通过！")

test_periodic_mapping()