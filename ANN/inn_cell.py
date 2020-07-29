from __future__ import print_function
import argparse
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.utils.data
import random
import numpy
import scipy.stats
#from torchvision import datasets, transforms
import matplotlib.pyplot as plt


from torch.utils.data.dataset import Dataset
from torch.utils.data import DataLoader



class MyDataset(Dataset):
    def __init__(self, X, y):
        self.data = X
        self.target = y
    
    def __getitem__(self, index):
        x = self.data[index]
        y = self.target[index]
        return x, y

    def __len__(self):
        return len(self.data)





class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()

        if (args.act == 'mixed'):
            self.layer1 = nn.Linear(NSUPRA, 5, bias=False)
            self.layer2 = nn.Linear(NDENDS-NSUPRA, 5, bias=False)
            self.soma_layer = nn.Linear(10, 1)
        else:
            self.layer1 = nn.Linear(NDENDS, 5, bias=False)
            self.soma_layer = nn.Linear(5, 1)

    @staticmethod
    def mysub(x):
           return torch.pow(x.clamp(0) + 2., 0.7) - 2.
           

    @staticmethod
    def mysupra(x):
           return 5*torch.sigmoid(x/5. - 5.)

    def forward(self, xinput):


        if (args.act == 'linear') :
            a = F.relu(self.layer1(xinput))
            c = a
        elif (args.act == 'sub') :
            a = Net.mysub(self.layer1(xinput))
            c = a
        elif (args.act == 'supra') :
            a = Net.mysupra(self.layer1(xinput))
            c = a
        else: # mixed
            in_supra = xinput[:, 0:NSUPRA];
            in_sub = xinput[:, NSUPRA:];
            a = Net.mysupra(self.layer1(in_supra))
            b = Net.mysub(self.layer2(in_sub))
            c = torch.cat((a,b), 1 ) ;

        soma_out = self.soma_layer(c)

        return soma_out


def train(args, model, device, train_loader, optimizer, epoch):
    model.train()

    tot_loss = 0

    for batch_idx, (data, target) in enumerate(train_loader):
        data, target = data.to(device), target.to(device)

        optimizer.zero_grad()
        output = model(data)

        loss = F.mse_loss(output, target, reduce=True)

        loss.backward()
        optimizer.step()

        for p in model.parameters():
            p.data.clamp_(0)
        
        tot_loss += loss.item()


        if (False):
            if batch_idx % args.log_interval == 0:
                print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                    epoch, batch_idx * len(data), len(train_loader.dataset),
                    100. * batch_idx / len(train_loader), loss.item()))

    tot_loss /= len(train_loader.dataset)
    #print('Epoch loss: %f'%(tot_loss))
    return tot_loss



def test(args, model, device, test_loader):
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for data, target in test_loader:
            data, target = data.to(device), target.to(device)
            output = model(data)
            #test_loss += F.mse_loss(output, target, reduction='sum').item() # sum up batch loss
            test_loss += F.mse_loss(output, target, reduce=True)
 
    test_loss /= len(test_loader.dataset)

    return test_loss



def pred(args, model, device, loader):
    model.eval()
    test_loss = 0

    # fill in with predicted values
    out_data = torch.zeros((len(loader.dataset),1))
    idx=0
    with torch.no_grad():
        for data, target in loader:
            data, target = data.to(device), target.to(device)
            out_data[idx:idx+len(data)] = model(data)
            #test_loss += F.mse_loss(output, target, reduction='sum').item() # sum up batch loss
            #test_loss += F.mse_loss(output, target, reduce=True)
            idx += len(data)

    #test_loss /= len(loader.dataset)
    # print("\nTest set: Average loss: %4f\n"%( test_loss) ) 
    return out_data



# Training settings
parser = argparse.ArgumentParser(description='PyTorch Example')
parser.add_argument('--epochs', type=int, default=20, metavar='N',
                    help='number of epochs to train (default: 10)')
parser.add_argument('--lr', type=float, default=0.01, metavar='LR',
                    help='learning rate (default: 0.01)')
parser.add_argument('--no-cuda', action='store_true', default=False,
                    help='disables CUDA training')
parser.add_argument('--seed', type=int, default=1, metavar='S',
                    help='random seed (default: 1)')
parser.add_argument('--log-interval', type=int, default=10, metavar='N',
                    help='how many batches to wait before logging training status')
parser.add_argument('--act', default='linear',
                    help='activation: linear, sub, supra')

parser.add_argument('--cell', default='1',
                    help='cell number to read data from')

parser.add_argument('--plot',  help='show plots?' , action='store_true')

parser.add_argument('--pfc',  help='Use PFC data  instead of HIPP ? ' , action='store_true')

args = parser.parse_args()


#### MAIN ######

torch.manual_seed(args.seed)
numpy.random.seed(args.seed)
random.seed(args.seed)


NBATCH = 20 


# NUMBER OF SUPRA dendrites for each of the 8 cells
SUPRA_NUMBERS = [ 106, 13, 10, 90 , 27, 34, 48, 43 ];


print("Cell number: %s | nSeed: %s | Activations: %s | Epochs: %d"%(args.cell, args.seed, args.act, args.epochs))
print("Plot graphs: %d"%(args.plot))


X =  numpy.loadtxt(open("./CSVDATA/x_data_%s.csv"%(args.cell), "rb"), delimiter=",", skiprows=0)
Y =  numpy.loadtxt(open("./CSVDATA/y_data_%s.csv"%(args.cell), "rb"), delimiter=",", skiprows=0)
NSUPRA= SUPRA_NUMBERS[int(args.cell)-1]
print("Reading xdata from  x_data_%s.csv with %d supra dendrites "%(args.cell, NSUPRA))


# number of samples
NSAMPLES = X.shape[0]
NSAMPLES -= NSAMPLES%NBATCH #trim data set to fit batch size

print("Total samples: %d"%(NSAMPLES))

X = X[0:NSAMPLES]
Y = Y[0:NSAMPLES]

NDENDS   = X.shape[1]

# ensure they are 32bit floats  or we ll get an error
all_x = torch.tensor(numpy.reshape(X, (NSAMPLES, NDENDS))).float()    
all_y = torch.tensor(numpy.reshape(Y, (NSAMPLES,1))).float()




#shuffle
p = numpy.random.permutation(NSAMPLES)


# Get approx 80% percent of data for the training set
# The size of input data must be a multiple of the batch size, so we trim some data:

NSPLIT = int(float(NSAMPLES)*0.80);
NSPLIT -= NSPLIT%NBATCH

x_train = all_x[p[0:NSPLIT]]
y_train = all_y[p[0:NSPLIT]]

print("Shape of training data:")
print(x_train.shape)


# Test set data
x_test = all_x[p[NSPLIT:NSAMPLES]]
y_test = all_y[p[NSPLIT:NSAMPLES]]

print("Shape of test data:")
print(x_test.shape)


use_cuda = False #not args.no_cuda and torch.cuda.is_available()

device = torch.device("cuda" if use_cuda else "cpu")

kwargs = {'num_workers': 1, 'pin_memory': True} if use_cuda else {}

train_dataset = MyDataset(x_train, y_train)
test_dataset = MyDataset(x_test, y_test)
pred_dataset = MyDataset(all_x[0:NSAMPLES], all_y)

train_loader = DataLoader(train_dataset, batch_size=NBATCH, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=NBATCH, shuffle=True)


model = Net().to(device)
#optimizer = optim.SGD(model.parameters(), lr=args.lr)
optimizer = optim.Adam(model.parameters(), lr=args.lr, betas=(0.9, 0.999))


error_log = numpy.zeros((args.epochs, 2))
for epoch in range(0, args.epochs):
    error_log[epoch, 0 ] = train(args, model, device, train_loader, optimizer, epoch)
    error_log[epoch, 1 ] = test(args, model, device, test_loader)
    if (epoch%10==9):
        print("[epoch %d] Train set  loss: %4f -  test: %f"%( epoch, error_log[epoch, 0 ], error_log[epoch, 1 ]) ) 

numpy.savetxt('./data/errors-%s-%d-%s.txt'%(args.cell, args.seed, args.act), error_log, fmt="%f")



#plot_dataset = pred_dataset;
#actual_data = all_y
x_data=x_test  
## make predictions from only the test data set
plot_dataset = test_dataset;
actual_data = y_test

plot_loader = DataLoader(plot_dataset, batch_size=NBATCH, shuffle=False)

predictions =  pred(args, model, device, plot_loader)

print('Saving predictions-%s-%d-%s.txt'%(args.cell, args.seed, args.act))
numpy.savetxt('./data/predictions-%s-%d-%s.txt'%(args.cell, args.seed, args.act), predictions, fmt="%f")

print('Saving actual-%s-%d.txt'%(args.cell, args.seed))
numpy.savetxt('./data/actual-%s-%d.txt'%(args.cell, args.seed), actual_data, fmt="%f")

print('Saving x-%s-%d.txt'%(args.cell, args.seed))
numpy.savetxt('./data/x-%s-%d.txt'%(args.cell, args.seed), x_data, fmt="%f")


print("Finito %s"%(args.act))


if (args.plot):
    plt.figure()

    plt.plot(error_log[:, 0])
    plt.plot(error_log[:, 1])
    plt.xlabel("Epoch")
    plt.ylabel("MSE Error")
    plt.legend(["Train error","Test error"] )

    plt.figure()
    plt.scatter( actual_data, predictions)
    pearsonr = scipy.stats.pearsonr(actual_data, predictions)
    plt.xlabel('Actual')
    plt.title("%s,  r=%f p=%E"%( args.act,  pearsonr[0], pearsonr[1]))
    print("%s,  r=%f p=%E"%( args.act,  pearsonr[0], pearsonr[1]))

    plt.show()

