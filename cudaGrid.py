# Import libraries for simulation
import tensorflow as tf
import numpy as np

# Imports for visualization
import PIL.Image


def DisplayArray(a, name, fmt='jpeg', rng=[0, 1]):
    """Display an array as a picture."""
    a = (a - rng[0])/float(rng[1] - rng[0])*255
    a = np.uint8(np.clip(a, 0, 255))
    with open("fig/{0}".format(name), "w") as f:
        PIL.Image.fromarray(a).save(f, "jpeg")


# sess = tf.Session()
sess = tf.InteractiveSession()

# Computational Convenience Functions


def make_kernel(a):
    """Transform a 2D array into a convolution kernel"""
    # https://www.youtube.com/watch?v=6yrPU8rYOhs
    a = np.asarray(a)
    a = a.reshape(list(a.shape) + [1, 1])
    return tf.constant(a, dtype=1)


def simple_conv(x, k):
    """A simplified 2D convolution operation"""
    x = tf.expand_dims(tf.expand_dims(x, 0), -1)
    y = tf.nn.depthwise_conv2d(x, k, [1, 1, 1, 1], padding='SAME')
    return y[0, :, :, 0]


# https://en.wikipedia.org/wiki/Discrete_Laplace_operator
def laplace(x):
    """Compute the 2D laplacian of an array"""
    laplace_k = make_kernel([[0.5, 1.0, 0.5],
                             [1.0, -6., 1.0],
                             [0.5, 1.0, 0.5]])
    return simple_conv(x, laplace_k)


# Define the PDE
# plate size, mm
w = h = 1000
# intervals in x-, y- directions, mm
dx = dy = 1
# Thermal diffusivity of steel, mm2.s-1
D = 0.1328
nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = dx2 * dy2 / (2 * D * (dx2 + dy2))
Tice, Tcore, Tmelt = np.float32(-20.0), np.float32(500.0), np.float32(0.0)
r, cx, cy = 16, w/2., h/2.
r2 = r**2

# Initial Conditions

# Set everything to Tice
u_init = Tice * np.ones((nx, ny), dtype=np.float32)
ut_init = Tice * np.ones((nx, ny), dtype=np.float32)

# Set root temp to Thot
for i in range(nx):
    for j in range(ny):
        p2 = (i*dx-cx)**2 + (j*dy-cy)**2
        if p2 < r2:
            u_init[i, j] = Tcore

DisplayArray(u_init, name="out-{0}.png".format(str(0).zfill(6), rng=[-0.1, 0.1], fmt='png'))

# Parameters:
# eps -- time resolution
# damping -- wave damping
eps = tf.placeholder(tf.float32, shape=())
damping = tf.placeholder(tf.float32, shape=())

# Create variables for simulation state
U = tf.Variable(u_init)
Ut = tf.Variable(ut_init)

# Discretized PDE update rules
U_ = U + eps * Ut
Ut_ = Ut + eps * (laplace(U) - damping * Ut)

# Operation to update the state
step = tf.group(
    U.assign(U_),
    Ut.assign(Ut_))

# Initialize state to initial conditions
tf.initialize_all_variables().run()

# Run 1000 steps of PDE
for i in range(1000):
    # Step simulation
    step.run({eps: 0.1, damping: 0.0})
    if (i % 10) == 0:
        DisplayArray(U.eval(), name="out-{0}.png".format(str(i).zfill(6), rng=[-0.1, 0.1], fmt='png'))
