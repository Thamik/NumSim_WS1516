from uq_sampling import SimFrame

if __name__ == '__main__':
        # only run when executed as a script
        frame = SimFrame(2000, 4)
        frame.generateSamples()
        frame.run()
