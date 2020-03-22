'''
This is definitely an overestimation, but i wanted to keep some track so i just quickly made this to have an idea.

There might be better ways of doing this when you have more of a finished product but this will do for now.

'''



filename = 'kththesis.tex'

actual_text = []
with open(filename) as f:
    lines = f.readlines()

    for i in lines:
        i = i.lstrip()
        if i.startswith("\\"):
            pass
        elif i.startswith('%'):
            pass
        else:
            i = i.strip('\n')

            actual_text.append(i)

actual_text = list(filter(None, actual_text))

count = 0
for i in actual_text:
    i = i.split(' ')
    count += len(i)

print(count)
