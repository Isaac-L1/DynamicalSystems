function send_complete(dest, task)

    mail = 'ielhome253@gmail.com';
    password = 'Bananas253';
    setpref('Internet','SMTP_Server','smtp.gmail.com');
    setpref('Internet','E_mail',mail);
    setpref('Internet','SMTP_Username',mail);
    setpref('Internet','SMTP_Password',password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.starttls.enable','true');
    %props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');
    sendmail(dest, 'Task Complete', ['Your task: ' task ' has completed.']);

end