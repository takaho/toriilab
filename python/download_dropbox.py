import sys
import dropbox
from dropbox.files import WriteMode
from dropbox.exceptions import ApiError, AuthError

AppKey = '51kyqv4dw7w28ck'
TOKEN = '5HWO0GV38y8AAAAAAAAMjw-fQr7Be7p8NHpvqi_toszRdm3U81KB-A0sgAoD6L58'
#TOKEN = 'Q4u8iEux9ZAAAAAAAAAAf7k-JWXvHDlRrKsVLbLmeKeydepkcCxui4RaACQjF9NS'

def list_folder(dbx, path):
    #with stopwatch('list_folder'):
    res = dbx.files_list_folder(path)
    print(res)
    entries = {}
    for entry in res.entries:
        entries[entry.name] = entry
        print(entry.name)
    return entries

if __name__ == '__main__':
    dbx = dropbox.Dropbox(TOKEN)

    # print(dbx)
    # try:
    #     dbx.users_get_current_account()
    # except:
    #     sys.exit("ERROR: Invalid access token")
    #     raise

    #list_folder(dbx, '/home/Fastq/fastq_181217_nextseq')
    list_folder(dbx, '/fastq_181217_nextseq') 
   #list_folder(dbx, '/fastq_181217_nextseq')
    # #md, res = dbx.files_download(path)
    # #data = res.content
    
