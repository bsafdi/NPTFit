while True:
    try:
        n = raw_input("Please enter an integer: ")
        n = int(n)
        #break
        pass
    except ValueError:
        print("No valid integer! Please try again ...")
print "Great, you successfully entered an integer!"
