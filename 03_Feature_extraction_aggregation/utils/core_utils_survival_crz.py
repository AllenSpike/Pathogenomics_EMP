import pdb
from argparse import Namespace
from collections import OrderedDict
import os
import pickle

from lifelines.utils import concordance_index
import numpy as np
from sksurv.metrics import concordance_index_censored

import torch

# NOTE: 这里使用model_clam
from datasets.dataset_generic import save_splits
from models.model_clam import CLAM_MB, CLAM_SB

from utils.utils_survival import *
from utils.loss_func import NLLSurvLoss

from utils.coattn_train_utils import *
from utils.cluster_train_utils import *

from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import auc as calc_auc
from sklearn.preprocessing import label_binarize

class Accuracy_Logger(object):
    """Accuracy logger"""

    def __init__(self, n_classes):
        super(Accuracy_Logger, self).__init__()
        self.n_classes = n_classes
        self.initialize()

    def initialize(self):
        self.data = [{"count": 0, "correct": 0} for i in range(self.n_classes)]

    def log(self, Y_hat, Y):
        Y_hat = int(Y_hat)
        Y = int(Y)
        self.data[Y]["count"] += 1
        self.data[Y]["correct"] += (Y_hat == Y)

    def log_batch(self, Y_hat, Y):
        Y_hat = np.array(Y_hat).astype(int)
        Y = np.array(Y).astype(int)
        for label_class in np.unique(Y):
            cls_mask = Y == label_class
            self.data[label_class]["count"] += cls_mask.sum()
            self.data[label_class]["correct"] += (Y_hat[cls_mask] == Y[cls_mask]).sum()

    def get_summary(self, c):
        count = self.data[c]["count"]
        correct = self.data[c]["correct"]

        if count == 0:
            acc = None
        else:
            acc = float(correct) / count

        return acc, correct, count

class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""

    def __init__(self, warmup=5, patience=15, stop_epoch=20, verbose=False):
        """
        Args:
            patience (int): How long to wait after last time validation loss improved.
                            Default: 20
            stop_epoch (int): Earliest epoch possible for stopping
            verbose (bool): If True, prints a message for each validation loss improvement.
                            Default: False
        """
        self.warmup = warmup
        self.patience = patience
        self.stop_epoch = stop_epoch
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf

    def __call__(self, epoch, val_loss, model, ckpt_name='checkpoint.pt'):

        score = -val_loss

        if epoch < self.warmup:
            pass
        elif self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model, ckpt_name)
        elif score < self.best_score:
            self.counter += 1
            print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience and epoch > self.stop_epoch:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model, ckpt_name)
            self.counter = 0

    def save_checkpoint(self, val_loss, model, ckpt_name):
        '''Saves model when validation loss decrease.'''
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')
        torch.save(model.state_dict(), ckpt_name)
        self.val_loss_min = val_loss


class Monitor_CIndex:
    """Early stops the training if validation loss doesn't improve after a given patience."""

    def __init__(self):
        """
        Args:
            patience (int): How long to wait after last time validation loss improved.
                            Default: 20
            stop_epoch (int): Earliest epoch possible for stopping
            verbose (bool): If True, prints a message for each validation loss improvement.
                            Default: False
        """
        self.best_score = None

    def __call__(self, val_cindex, model, ckpt_name: str = 'checkpoint.pt'):

        score = val_cindex

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(model, ckpt_name)
        elif score > self.best_score:
            self.best_score = score
            self.save_checkpoint(model, ckpt_name)
        else:
            pass

    def save_checkpoint(self, model, ckpt_name):
        '''Saves model when validation loss decrease.'''
        torch.save(model.state_dict(), ckpt_name)


# SECTION：模型训练
def train_survival(datasets: tuple, cur: int, args: Namespace):
    """
        train for a single fold
    """
    print('\nTraining Fold {}!'.format(cur))
    writer_dir = os.path.join(args.results_dir, str(cur))
    if not os.path.isdir(writer_dir):
        os.mkdir(writer_dir)

    if args.log_data:
        from tensorboardX import SummaryWriter
        writer = SummaryWriter(writer_dir, flush_secs=15)

    else:
        writer = None

    print('\nInit train/val/test splits...', end=' ')
    train_split, val_split, test_split = datasets
    save_splits(datasets, ['train', 'val', 'test'], os.path.join(args.results_dir, 'splits_{}.csv'.format(cur)))
    print('Done!')
    print("Training on {} samples".format(len(train_split)))
    print("Validating on {} samples".format(len(val_split)))
    print("Testing on {} samples".format(len(test_split)))

    print('\nInit loss function...', end=' ')
    if args.task == 'survival':
        if args.bag_loss == 'nll_surv': # NOTE:bag loss
            loss_fn = NLLSurvLoss(alpha=args.alpha_surv)
        elif args.bag_loss == 'ce_surv':
            loss_fn = CrossEntropySurvLoss(alpha=args.alpha_surv)
        elif args.bag_loss == 'ce':
            loss_fn = nn.CrossEntropyLoss()
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError
    print('Done!')

    print('\nInit Model...', end=' ')

    ### CLAM
    if args.model_type in ['clam_sb', 'clam_mb']:

        model_dict = {"dropout": args.drop_out, 'n_classes': args.n_classes}

        if args.subtyping:
            model_dict.update({'subtyping': True})

        if args.B > 0:
            model_dict.update({'k_sample': args.B})

        if args.inst_loss == 'svm': # NOTE: instance loss
            from topk.svm import SmoothTop1SVM
            instance_loss_fn = SmoothTop1SVM(n_classes=2)
            if device.type == 'cuda':
                instance_loss_fn = instance_loss_fn.cuda()
        else:
            instance_loss_fn = nn.CrossEntropyLoss()

        if args.model_type == 'clam_sb':
            model = CLAM_SB(**model_dict, instance_loss_fn=instance_loss_fn)
        elif args.model_type == 'clam_mb':
            model = CLAM_MB(**model_dict, instance_loss_fn=instance_loss_fn)
        else:
            raise NotImplementedError

    if hasattr(model, "relocate"):
        model.relocate()
    else:
        model = model.to(torch.device('cuda'))
    print('Done!')
    print_network(model)

    print('\nInit optimizer ...', end=' ')
    optimizer = get_optim(model, args)
    print('Done!')

    print('\nInit Loaders...', end=' ')
    train_loader = get_split_loader_surv(split_dataset=train_split, training=True, testing=args.testing,
                                         weighted=args.weighted_sample, mode=args.mode, batch_size=1)
    val_loader = get_split_loader_surv(split_dataset=val_split, testing=args.testing,
                                       mode=args.mode, batch_size=1)
    test_loader = get_split_loader_surv(split_dataset=test_split, testing=args.testing,
                                        mode=args.mode, batch_size=1)
    print('Done!')

    print('\nSetup EarlyStopping...', end=' ')
    if args.early_stopping:
        early_stopping = EarlyStopping(warmup=0, patience=10, stop_epoch=20, verbose=True)
    else:
        early_stopping = None

    print('\nSetup Validation C-Index Monitor...', end=' ')
    monitor_cindex = Monitor_CIndex()
    print('Done!')

    for epoch in range(args.max_epochs):
        if args.task == 'survival':
            if args.model_type in ['clam_sb', 'clam_mb'] and not args.no_inst_cluster:
                train_loop_clam(epoch, model, train_loader, optimizer, args.n_classes, args.bag_weight, writer, loss_fn)
                stop = validate_clam(cur, epoch, model, val_loader, args.n_classes,
                                     early_stopping, writer, loss_fn, args.results_dir)
            if stop:
                break

    if args.early_stopping:
        model.load_state_dict(torch.load(os.path.join(args.results_dir, "s_{}_checkpoint.pt".format(cur)),
                                         weights_only=True))  # NOTE:官方警告
    else:
        torch.save(model.state_dict(), os.path.join(args.results_dir, "s_{}_checkpoint.pt".format(cur)))

    results_val_dict, val_cindex, val_error, val_auc, _ = summary_survival(model, val_loader, args.n_classes)
    print('Val c-Index: {:.4f}, Val error: {:.4f}, ROC AUC: {:.4f}'.format(val_cindex, val_error, val_auc)) # NOTE:验证集C-index
    results_test_dict, test_cindex, test_error, test_auc, _ = summary_survival(model, test_loader, args.n_classes) # NOTE: 测试集C-index
    print('Test c-Index: {:.4f}, Test error: {:.4f}, Test AUC: {:.4f}'.format(test_cindex, test_error, test_auc))

    writer.close()

    return results_test_dict, test_cindex, val_cindex, test_auc, val_auc # NOTE:返回测试集结果、测试集C-index、验证集C-index、测试集AUC、验证集AUC

def summary_survival(model, loader, n_classes):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    acc_logger = Accuracy_Logger(n_classes=n_classes)
    model.eval()
    test_loss = 0.
    test_error = 0.

    all_probs = np.zeros((len(loader), n_classes))
    all_labels = np.zeros(len(loader))

    all_risk_scores = np.zeros((len(loader)))
    all_censorships = np.zeros((len(loader)))
    all_event_times = np.zeros((len(loader)))

    slide_ids = loader.dataset.slide_data['slide_id']
    patient_results = {}

    for batch_idx, (data_WSI, data_omic, y_disc, event, c) in enumerate(loader):

        data = data_WSI.to(device)
        label = y_disc.to(device)
        event_time = event.to(device)
        censor = c.to(device)

        slide_id = slide_ids.iloc[batch_idx]

        with torch.no_grad():
            logits, Y_prob, Y_hat, _, _ = model(data)

        acc_logger.log(Y_hat, label)
        probs = Y_prob.cpu().numpy()
        all_probs[batch_idx] = probs
        all_labels[batch_idx] = label.item()

        hazards = torch.sigmoid(logits)
        survival = torch.cumprod(1 - hazards, dim=1)
        risk = -torch.sum(survival, dim=1).detach().cpu().numpy()

        # pdb.set_trace()

        # event_time = np.asscalar(event_time)
        # censor = np.asscalar(censor)
        event_time = event_time.item()
        censor = censor.item()  # NOTE：np.asscalar已弃用

        all_risk_scores[batch_idx] = risk
        all_censorships[batch_idx] = censor
        all_event_times[batch_idx] = event_time
        patient_results.update({slide_id: {'slide_id': np.array(slide_id), 'prob': probs, 'risk': risk, 'disc_label': y_disc.item(),
                                           'survival': event_time, 'censorship': censor}})

    c_index = \
    concordance_index_censored((1 - all_censorships).astype(bool), all_event_times, all_risk_scores, tied_tol=1e-08)[0]

    test_error /= len(loader)

    if n_classes == 2:
        auc = roc_auc_score(all_labels, all_probs[:, 1])
    else:
        aucs = []
        binary_labels = label_binarize(all_labels, classes=[i for i in range(n_classes)])
        for class_idx in range(n_classes):
            if class_idx in all_labels:
                fpr, tpr, _ = roc_curve(binary_labels[:, class_idx], all_probs[:, class_idx])
                aucs.append(calc_auc(fpr, tpr))
            else:
                aucs.append(float('nan'))

        auc = np.nanmean(np.array(aucs))

    return patient_results, c_index, test_error, auc, acc_logger


def train_loop_clam(epoch, model, loader, optimizer, n_classes, bag_weight, writer=None, loss_fn=None):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model.train()

    acc_logger = Accuracy_Logger(n_classes=n_classes)
    inst_logger = Accuracy_Logger(n_classes=n_classes)

    train_loss = 0.
    train_error = 0.
    train_inst_loss = 0.
    inst_count = 0

    all_risk_scores = []
    all_censorships = []
    all_event_times = []

    print('\n')

    for batch_idx, (data_WSI, data_omic, y_disc, event, c) in enumerate(loader):

        data = data_WSI.to(device)
        label = y_disc.to(device)
        event_time = event.to(device)
        censor = c.to(device)

        logits, Y_prob, Y_hat, _, instance_dict = model(data, label=label, instance_eval=True, n_classes=n_classes)
        acc_logger.log(Y_hat, label)

        hazards = torch.sigmoid(logits)
        survival = torch.cumprod(1 - hazards, dim=1)
        risk = -torch.sum(survival, dim=1).detach().cpu().numpy()

        if isinstance(loss_fn, NLLSurvLoss):
            loss = loss_fn(h=logits, y=label, t=event_time, c=censor)
        elif isinstance(loss_fn, CrossEntropySurvLoss):
            loss = loss_fn(hazards=hazards, Y=label, c=censor, S=None)
        elif isinstance(loss_fn, nn.CrossEntropyLoss):
            loss = loss_fn(input=logits, target=label.long())
        loss_value = loss.item() # bag loss

        instance_loss = instance_dict['instance_loss']
        inst_count += 1
        instance_loss_value = instance_loss.item() # instance loss

        train_inst_loss += instance_loss_value
        train_loss += loss_value
        total_loss = bag_weight * loss + (1 - bag_weight) * instance_loss  # total loss

        inst_preds = instance_dict['inst_preds']
        inst_labels = instance_dict['inst_labels']
        inst_logger.log_batch(inst_preds, inst_labels)

        if (batch_idx + 1) % 20 == 0:
            print('batch {}, loss: {:.4f}, instance_loss: {:.4f}, weighted_loss: {:.4f}, '.format(batch_idx, loss_value,
                                                                                                      instance_loss_value,
                                                                                                      total_loss.item()) +
                      'label: {}, event_time: {:.4f}, risk: {:.4f}, bag_size: {}'.format(label.item(), float(event_time.detach().cpu()[0]), float(risk[0]), data.size(0)))

        error = calculate_error(Y_hat, label)
        train_error += error

        # backward pass
        total_loss.backward()

        # step
        optimizer.step()
        optimizer.zero_grad()

        all_risk_scores.append(risk)
        all_censorships.append(censor.detach().cpu().numpy())
        all_event_times.append(event_time.detach().cpu().numpy())

    # calculate loss and error for epoch
    train_loss /= len(loader)
    train_error /= len(loader)

    if inst_count > 0:
        train_inst_loss /= inst_count
        print('\n')
        for i in range(2):
            acc, correct, count = inst_logger.get_summary(i)
            print('class {} clustering acc {}: correct {}/{}'.format(i, acc, correct, count))

    all_risk_scores = np.concatenate(all_risk_scores)
    all_censorships = np.concatenate(all_censorships)
    all_event_times = np.concatenate(all_event_times)
    c_index = concordance_index_censored((1 - all_censorships).astype(bool), all_event_times, all_risk_scores, tied_tol=1e-08)[0]

    print('Epoch: {}, train_loss: {:.4f}, train_clustering_loss:  {:.4f}, train_error: {:.4f}, train_c_index: {:.4f}'.format(epoch, train_loss, train_inst_loss, train_error, c_index))

    for i in range(n_classes):
        acc, correct, count = acc_logger.get_summary(i)
        print('class {}: acc {}, correct {}/{}'.format(i, acc, correct, count))
        if writer and acc is not None:
            writer.add_scalar('train/class_{}_acc'.format(i), acc, epoch)

    if writer:
        writer.add_scalar('train/loss', train_loss, epoch)
        writer.add_scalar('train/error', train_error, epoch)
        writer.add_scalar('train/clustering_loss', train_inst_loss, epoch)


def validate_clam(cur, epoch, model, loader, n_classes, early_stopping=None, writer=None, loss_fn=None,
                  results_dir=None):

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model.eval()

    acc_logger = Accuracy_Logger(n_classes=n_classes)
    inst_logger = Accuracy_Logger(n_classes=n_classes)
    val_loss = 0.
    val_error = 0.

    val_inst_loss = 0.
    val_inst_acc = 0.
    inst_count = 0

    all_risk_scores = np.zeros((len(loader)))
    all_censorships = np.zeros((len(loader)))
    all_event_times = np.zeros((len(loader)))

    prob = np.zeros((len(loader), n_classes))
    labels = np.zeros(len(loader))
    sample_size = model.k_sample

    with torch.no_grad():

        for batch_idx, (data_WSI, data_omic, y_disc, event, c) in enumerate(loader):

            data = data_WSI.to(device)
            label = y_disc.to(device)
            event_time = event.to(device)
            censor = c.to(device)

            logits, Y_prob, Y_hat, _, instance_dict = model(data, label=label, instance_eval=True, n_classes=n_classes)
            acc_logger.log(Y_hat, label)

            hazards = torch.sigmoid(logits)
            survival = torch.cumprod(1 - hazards, dim=1)
            risk = -torch.sum(survival, dim=1).detach().cpu().numpy()

            if isinstance(loss_fn, NLLSurvLoss):
                loss = loss_fn(h=logits, y=label, t=event_time, c=censor)
            elif isinstance(loss_fn, CrossEntropySurvLoss):
                loss = loss_fn(hazards=hazards, Y=label, c=censor,S=None)
            elif isinstance(loss_fn, nn.CrossEntropyLoss):
                loss = loss_fn(logits, label.long())
            val_loss = loss.item()  # bag loss

            instance_loss = instance_dict['instance_loss']
            inst_count += 1
            instance_loss_value = instance_loss.item() # instance loss

            val_inst_loss += instance_loss_value
            inst_preds = instance_dict['inst_preds']
            inst_labels = instance_dict['inst_labels']
            inst_logger.log_batch(inst_preds, inst_labels)

            prob[batch_idx] = Y_prob.cpu().numpy()
            labels[batch_idx] = label.item()

            error = calculate_error(Y_hat, label)
            val_error += error

            all_risk_scores[batch_idx] = risk
            all_censorships[batch_idx] = censor.detach().cpu().numpy()
            all_event_times[batch_idx] = event_time.detach().cpu().numpy()

    val_error /= len(loader)
    val_loss /= len(loader)
    c_index = concordance_index_censored((1 - all_censorships).astype(bool), all_event_times, all_risk_scores, tied_tol=1e-08)[0]

    if n_classes == 2:
        auc = roc_auc_score(labels, prob[:, 1])
        aucs = []
    else:
        aucs = []
        binary_labels = label_binarize(labels, classes=[i for i in range(n_classes)])
        for class_idx in range(n_classes):
            if class_idx in labels:
                fpr, tpr, _ = roc_curve(binary_labels[:, class_idx], prob[:, class_idx])
                aucs.append(calc_auc(fpr, tpr))
            else:
                aucs.append(float('nan'))

        auc = np.nanmean(np.array(aucs))

    print('\nVal Set, val_loss: {:.4f}, val_error: {:.4f}, c_index, {:.4f}, auc: {:.4f}'.format(val_loss, val_error, c_index, auc))

    if inst_count > 0:
        val_inst_loss /= inst_count
        for i in range(2):
            acc, correct, count = inst_logger.get_summary(i)
            print('class {} clustering acc {}: correct {}/{}'.format(i, acc, correct, count))

    if writer:
        writer.add_scalar('val/loss', val_loss, epoch)
        writer.add_scalar('val/auc', auc, epoch)
        writer.add_scalar('val/error', val_error, epoch)
        writer.add_scalar('val/inst_loss', val_inst_loss, epoch)
        writer.add_scalar('val/c-index', c_index, epoch)

    for i in range(n_classes):
        acc, correct, count = acc_logger.get_summary(i)
        print('class {}: acc {}, correct {}/{}'.format(i, acc, correct, count))

        if writer and acc is not None:
            writer.add_scalar('val/class_{}_acc'.format(i), acc, epoch)

    if early_stopping:
        assert results_dir
        early_stopping(epoch, val_loss, model, ckpt_name=os.path.join(results_dir, "s_{}_checkpoint.pt".format(cur)))

        if early_stopping.early_stop:
            print("Early stopping")
            return True

    return False